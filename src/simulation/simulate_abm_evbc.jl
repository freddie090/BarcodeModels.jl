function _lineage_df(records::Vector{LineageRecord}, rep::Int64)
    if isempty(records)
        return DataFrame(id = Int64[], parent_id = Int64[], birth_time = Float64[], rep = Int64[])
    end
    ids = [r.id for r in records]
    parent_ids = [r.parent_id for r in records]
    birth_times = [r.birth_time for r in records]
    return DataFrame(id = ids, parent_id = parent_ids, birth_time = birth_times, rep = fill(rep, length(records)))
end

function _run_abm_passage_experiment!(
    model::ResPop_ABM_EvBC,
    cells::Vector{CancerCellEvBC};
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
    next_cell_id::Int64,
    lineage_records::Vector{LineageRecord},
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
        state = ResPopABMEvBCState(cells, next_cell_id, lineage_records)
        kmc_out = run_model_core_abm(model, state, sim; treat = treat)
        next_cell_id = state.next_cell_id
        lineage_records = state.lineage_records

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
                    extend_with_dead_cells!(cells, Nbuff, make_dead_cell_evbc)
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
                                         samp_cell_lin_df_vec = samp_cell_lin_df_vec)
                    live_cells = alive_cells(cells)
                    cells = sample(live_cells, Nseed, replace = false)
                    extend_with_dead_cells!(cells, Nbuff, make_dead_cell_evbc)
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
        "lineage_records" => lineage_records,
        "next_cell_id" => next_cell_id,
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

function _run_abm_passage_experiment!(
    model::ResDmg_ABM_EvBC,
    cells::Vector{ResDmgCellEvBC};
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
    next_cell_id::Int64,
    lineage_records::Vector{LineageRecord},
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
        state = ResDmgABMEvBCState(cells, next_cell_id, lineage_records)
        kmc_out = run_model_core_abm(model, state, sim; treat = treat)
        next_cell_id = state.next_cell_id
        lineage_records = state.lineage_records

        kmc_last_t = last(kmc_out.tvec)
        isfinite(kmc_last_t) || error("ABM produced a non-finite time value.")

        curr_t_candidate = min(max(kmc_last_t, curr_t), next_t)

        live_count = length(alive_cells(cells))
        if last(kmc_out.Nvec) >= Nmax
            if last(kmc_out.Pvec) >= n_pass_eff
                _record_resdmg_abm_outputs!(kmc_out, cells, rep, curr_P, cell_lin_df_vec,
                                            Nvec, nS_vec, nDS_vec, nDR_vec, nR_vec, tvec, Pvec,
                                            sub_sample_cells = sub_sample_cells, K = K,
                                            samp_cell_lin_df_vec = samp_cell_lin_df_vec)
                curr_t = min(max(kmc_last_t, curr_t), tmax)
                break
            else
                _record_resdmg_abm_outputs!(kmc_out, cells, rep, curr_P, cell_lin_df_vec,
                                            Nvec, nS_vec, nDS_vec, nDR_vec, nR_vec, tvec, Pvec,
                                            sub_sample_cells = sub_sample_cells, K = K,
                                            samp_cell_lin_df_vec = samp_cell_lin_df_vec)
                if live_count < Nseed
                    break
                else
                    live_cells = alive_cells(cells)
                    cells = sample(live_cells, Nseed, replace = false)
                    extend_with_dead_cells!(cells, Nbuff, make_dead_resdmg_cell_evbc)
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
                                        samp_cell_lin_df_vec = samp_cell_lin_df_vec)
            curr_t = min(max(kmc_last_t, curr_t), tmax)
            break
        elseif kmc_last_t >= (tmax - progress_tol)
            _record_resdmg_abm_outputs!(kmc_out, cells, rep, curr_P, cell_lin_df_vec,
                                        Nvec, nS_vec, nDS_vec, nDR_vec, nR_vec, tvec, Pvec,
                                        sub_sample_cells = sub_sample_cells, K = K,
                                        samp_cell_lin_df_vec = samp_cell_lin_df_vec)
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
                                                samp_cell_lin_df_vec = samp_cell_lin_df_vec)
                    live_cells = alive_cells(cells)
                    cells = sample(live_cells, Nseed, replace = false)
                    extend_with_dead_cells!(cells, Nbuff, make_dead_resdmg_cell_evbc)
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
        "lineage_records" => lineage_records,
        "next_cell_id" => next_cell_id,
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

function _simulate_experiment_abm(model::ResPop_ABM_EvBC, exp::ExperimentParams; kwargs...)
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
    dt_save_at = _kw(kwargs, :dt_save_at, model.abm.dt_save_at)
    run_IC = _kw(kwargs, :run_IC, exp.run_IC)
    run_colony = _kw(kwargs, :run_colony, exp.run_colony)

    @assert !(run_IC || run_colony) "run_IC and run_colony are not implemented yet for ResPop_ABM_EvBC."

    _validate_tmax_vector_constraints(exp.tmax, exp.t_Pass)
    _validate_tmax_length(exp.tmax, n_rep)

    n_pass_eff = exp.tmax isa AbstractVector ? 1 : (length(_passage_times(exp.t_Pass, Float64(exp.tmax))) + 1)

    model_eff = _with_drug_effect(model, de)
    base_model_eff = ResPop_ABM(model_eff.params; abm = model_eff.abm)
    rep_cells_base = _expand_split_cells_abm(base_model_eff, exp, n_rep;
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
    lineage_df_outs = DataFrame[]
    sim_dfs = DataFrame[]

    nseed_last = _nseed_last(exp.Nseed)
    for i in 1:n_rep
        rep_tmax = _replicate_tmax(exp.tmax, n_rep, i)
        rep_cells_evbc, next_cell_id, lineage_records = _to_evbc_cells(rep_cells_base[i]; t0 = 0.0)
        extend_with_dead_cells!(rep_cells_evbc, model.abm.Nbuff, make_dead_cell_evbc)

        sim = _run_abm_passage_experiment!(
            model_eff, rep_cells_evbc;
            t0 = 0.0, tmax = rep_tmax, t_Pass = exp.t_Pass,
            Nseed = nseed_last, Nmax = exp.Nmax, Cc = exp.Cc,
            treat_ons = exp.treat_ons, treat_offs = exp.treat_offs,
            dt_save_at = dt_save_at, Nbuff = model.abm.Nbuff,
            R_real = R_real, t_frac = t_frac, rep = i,
            treat = drug_treatment, drug_effect = de,
            next_cell_id = next_cell_id,
            lineage_records = lineage_records,
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
        push!(lineage_df_outs, _lineage_df(sim["lineage_records"], i))
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
                    t_outs[idx] = idx > 1 ? t_outs[idx - 1] : 0.0
                    u_outs[idx] = idx > 1 ? u_outs[idx - 1] : 0.0
                end
            end

            append!(fin_t_outs, t_outs)
            append!(fin_u_outs, u_outs)
        end
    end

    fin_t_outs = round.(fin_t_outs; digits = 0)
    fin_lin_df = join_dfs(lin_df_outs, "bc")
    fin_sol_df = isempty(sim_dfs) ? DataFrame() : vcat(sim_dfs...)
    fin_lineage_df = isempty(lineage_df_outs) ? DataFrame() : vcat(lineage_df_outs...)

    out = Dict{String, Any}(
        "lin_df" => fin_lin_df,
        "sol_df" => fin_sol_df,
        "lineage_df" => fin_lineage_df
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

function _simulate_experiment_abm(model::ResDmg_ABM_EvBC, exp::ExperimentParams; kwargs...)
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
    dt_save_at = _kw(kwargs, :dt_save_at, model.abm.dt_save_at)
    run_IC = _kw(kwargs, :run_IC, exp.run_IC)
    run_colony = _kw(kwargs, :run_colony, exp.run_colony)

    @assert !(run_IC || run_colony) "run_IC and run_colony are not implemented yet for ResDmg_ABM_EvBC."

    _validate_tmax_vector_constraints(exp.tmax, exp.t_Pass)
    _validate_tmax_length(exp.tmax, n_rep)

    n_pass_eff = exp.tmax isa AbstractVector ? 1 : (length(_passage_times(exp.t_Pass, Float64(exp.tmax))) + 1)

    model_eff = _with_drug_effect(model, de)
    base_model_eff = ResDmg_ABM(model_eff.params; abm = model_eff.abm)
    rep_cells_base = _expand_split_cells_abm(base_model_eff, exp, n_rep;
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
    lineage_df_outs = DataFrame[]
    sim_dfs = DataFrame[]

    nseed_last = _nseed_last(exp.Nseed)
    for i in 1:n_rep
        rep_tmax = _replicate_tmax(exp.tmax, n_rep, i)
        rep_cells_evbc, next_cell_id, lineage_records = _to_evbc_cells(rep_cells_base[i]; t0 = 0.0)
        extend_with_dead_cells!(rep_cells_evbc, model.abm.Nbuff, make_dead_resdmg_cell_evbc)

        sim = _run_abm_passage_experiment!(
            model_eff, rep_cells_evbc;
            t0 = 0.0, tmax = rep_tmax, t_Pass = exp.t_Pass,
            Nseed = nseed_last, Nmax = exp.Nmax, Cc = exp.Cc,
            treat_ons = exp.treat_ons, treat_offs = exp.treat_offs,
            dt_save_at = dt_save_at, Nbuff = model.abm.Nbuff,
            R_real = R_real, t_frac = t_frac, rep = i,
            treat = drug_treatment, drug_effect = de,
            next_cell_id = next_cell_id,
            lineage_records = lineage_records,
            sub_sample_cells = sub_sample_cells, K = K
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
        push!(lineage_df_outs, _lineage_df(sim["lineage_records"], i))
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
                    t_outs[idx] = idx > 1 ? t_outs[idx - 1] : 0.0
                    u_outs[idx] = idx > 1 ? u_outs[idx - 1] : 0.0
                end
            end

            append!(fin_t_outs, t_outs)
            append!(fin_u_outs, u_outs)
        end
    end

    fin_t_outs = round.(fin_t_outs; digits = 0)
    fin_lin_df = join_dfs(lin_df_outs, "bc")
    fin_sol_df = isempty(sim_dfs) ? DataFrame() : vcat(sim_dfs...)
    fin_lineage_df = isempty(lineage_df_outs) ? DataFrame() : vcat(lineage_df_outs...)

    out = Dict{String, Any}(
        "lin_df" => fin_lin_df,
        "sol_df" => fin_sol_df,
        "lineage_df" => fin_lineage_df
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
