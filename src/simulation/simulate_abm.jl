function _record_abm_outputs!(kmc_out, cells, rep::Int64, curr_P::Int64,
    cell_lin_df_vec::Vector{DataFrame},
    Nvec::Vector{Int64}, nS_vec::Vector{Int64}, nR_vec::Vector{Int64}, nE_vec::Vector{Int64},
    tvec::Vector{Float64}, Pvec::Vector{Int64};
    sub_sample_cells::Bool = false, K::Int64 = 0,
    samp_cell_lin_df_vec::Vector{DataFrame} = DataFrame[],
    cond::String = "DT")

    live_cells = alive_cells(cells)
    bc_df = get_counts(live_cells, string(cond, rep, "_P", curr_P))
    push!(cell_lin_df_vec, bc_df)

    if sub_sample_cells
        if K <= length(live_cells)
            samp_cells = sample(live_cells, K, replace = false)
            push!(samp_cell_lin_df_vec, get_counts(samp_cells, string(cond, rep, "_P", curr_P)))
        else
            push!(samp_cell_lin_df_vec, get_counts(live_cells, string(cond, rep, "_P", curr_P)))
        end
    end

    update_track_vec!(kmc_out, Nvec, nS_vec, nR_vec, nE_vec, tvec, Pvec)
end

function _append_invivo_abm_eg_outputs!(kmc_out,
    nS_EG0_vec::Vector{Int64}, nS_EG1_vec::Vector{Int64},
    nR_EG0_vec::Vector{Int64}, nR_EG1_vec::Vector{Int64},
    nE_EG0_vec::Vector{Int64}, nE_EG1_vec::Vector{Int64})

    append!(nS_EG0_vec, kmc_out.S_EG0_vec)
    append!(nS_EG1_vec, kmc_out.S_EG1_vec)
    append!(nR_EG0_vec, kmc_out.R_EG0_vec)
    append!(nR_EG1_vec, kmc_out.R_EG1_vec)
    append!(nE_EG0_vec, kmc_out.E_EG0_vec)
    append!(nE_EG1_vec, kmc_out.E_EG1_vec)
    nothing
end

function _append_invivo_abm_solution_outputs!(kmc_out,
    Nvec::Vector{Int64}, nS_vec::Vector{Int64}, nR_vec::Vector{Int64}, nE_vec::Vector{Int64},
    tvec::Vector{Float64}, Pvec::Vector{Int64},
    nS_EG0_vec::Vector{Int64}, nS_EG1_vec::Vector{Int64},
    nR_EG0_vec::Vector{Int64}, nR_EG1_vec::Vector{Int64},
    nE_EG0_vec::Vector{Int64}, nE_EG1_vec::Vector{Int64};
    t_offset::Float64, passage::Int64)

    append!(Nvec, kmc_out.Nvec)
    append!(nS_vec, kmc_out.Svec)
    append!(nR_vec, kmc_out.Rvec)
    append!(nE_vec, kmc_out.Evec)
    append!(tvec, kmc_out.tvec .+ t_offset)
    append!(Pvec, fill(passage, length(kmc_out.Pvec)))
    _append_invivo_abm_eg_outputs!(kmc_out,
                                   nS_EG0_vec, nS_EG1_vec,
                                   nR_EG0_vec, nR_EG1_vec,
                                   nE_EG0_vec, nE_EG1_vec)
    nothing
end

function _invivo_abm_sol_df(sim::Dict; cond::String, rep::Int64)
    return DataFrame(
        t = sim["tvec"],
        nS_EG0 = sim["nS_EG0_vec"],
        nS_EG1 = sim["nS_EG1_vec"],
        nR_EG0 = sim["nR_EG0_vec"],
        nR_EG1 = sim["nR_EG1_vec"],
        nE_EG0 = sim["nE_EG0_vec"],
        nE_EG1 = sim["nE_EG1_vec"],
        nS = sim["nS_vec"],
        nR = sim["nR_vec"],
        nE = sim["nE_vec"],
        n_EG0 = sim["nS_EG0_vec"] .+ sim["nR_EG0_vec"] .+ sim["nE_EG0_vec"],
        n_EG1 = sim["nS_EG1_vec"] .+ sim["nR_EG1_vec"] .+ sim["nE_EG1_vec"],
        N = sim["Nvec"],
        cond = cond,
        rep = rep,
        passage = sim["Pvec"]
    )
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

_unique_live_barcodes(cells) = unique([cell.barcode for cell in cells if cell.alive])

function _sample_uniform_barcodes(bcs::Vector{Float64}, n_target::Int64)
    n_target_clamped = clamp(n_target, 0, length(bcs))
    n_target_clamped == 0 && return Float64[]
    return sample(bcs, n_target_clamped, replace = false)
end

function _assign_resistance_by_barcode!(cells::Vector{CancerCell}, rho::Float64, n0::Int64)
    bcs = _unique_live_barcodes(cells)
    n_target = Int64(round(rho * n0))
    selected_bcs = Set(_sample_uniform_barcodes(bcs, n_target))
    for cell in cells
        cell.alive || continue
        cell.R = cell.barcode in selected_bcs
    end
    return nothing
end

function _assign_resistance_by_barcode!(cells::Vector{ResDmgCell}, rho::Float64, n0::Int64)
    bcs = _unique_live_barcodes(cells)
    n_target = Int64(round(rho * n0))
    selected_bcs = Set(_sample_uniform_barcodes(bcs, n_target))
    for cell in cells
        cell.alive || continue
        cell.R = cell.barcode in selected_bcs
    end
    return nothing
end

function _assign_resistance_by_barcode!(cells::Vector{InVivoCancerCell}, rho::Float64, n0::Int64)
    bcs = _unique_live_barcodes(cells)
    n_target = Int64(round(rho * n0))
    selected_bcs = Set(_sample_uniform_barcodes(bcs, n_target))
    for cell in cells
        cell.alive || continue
        cell.R = cell.barcode in selected_bcs
    end
    return nothing
end

function _assign_engraftment_by_barcode!(cells::Vector{InVivoCancerCell}, fEG1::Float64, n0::Int64)
    bcs = _unique_live_barcodes(cells)
    n_target = Int64(round(fEG1 * n0))
    selected_bcs = Set(_sample_uniform_barcodes(bcs, n_target))
    for cell in cells
        cell.alive || continue
        cell.EG = cell.barcode in selected_bcs
    end
    return nothing
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
    R_real::String = "b",
    drug_effect::Symbol = model.params.drug_effect,
    skew_lib::Bool = model.abm.skew_lib,
    use_lib_probs::Bool = model.abm.use_lib_probs,
    split_after_barcoding::Bool = model.abm.split_after_barcoding,
    bc_unif::Float64 = model.abm.bc_unif,
    Nbc::Int64 = model.abm.Nbc,
    bc_probs::Vector{Float64} = model.abm.bc_probs,
    dt_save_at::Float64 = model.abm.dt_save_at,
    t_frac::Float64 = model.abm.t_frac)

    Nbuff = model.abm.Nbuff
    barcode_kwargs = (; skew_lib = skew_lib, use_lib_probs = use_lib_probs,
                       bc_unif = bc_unif, Nbc = Nbc, bc_probs = bc_probs)

    if split_after_barcoding
        final_seed = _nseed_last(exp.Nseed)
        n_rep * final_seed <= exp.n0 || error("split_after_barcoding requires n0 >= n_rep*Nseed (n0=$(exp.n0), n_rep=$(n_rep), Nseed=$(final_seed)).")

        exp_cells = seed_cells(exp.n0, 0.0, Nbuff;
                               barcode_kwargs...)
        _assign_resistance_by_barcode!(exp_cells, model.params.rho, exp.n0)

        rep_cells = sample(alive_cells(exp_cells), n_rep * final_seed, replace = false)
        rep_cells = reshape(rep_cells, (final_seed, n_rep))

        fin_rep_cells = Vector{Vector{CancerCell}}(undef, n_rep)
        for i in 1:n_rep
            fin_rep_cells[i] = collect(rep_cells[:, i])
            extend_with_dead_cells!(fin_rep_cells[i], Nbuff, make_dead_cell)
        end
        return fin_rep_cells
    end

    exp_cells = seed_cells(exp.n0, model.params.rho, Nbuff;
                           barcode_kwargs...)

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
    use_lib_probs = _kw(kwargs, :use_lib_probs, model.abm.use_lib_probs)
    split_after_barcoding = _kw(kwargs, :split_after_barcoding, model.abm.split_after_barcoding)
    bc_unif = _kw(kwargs, :bc_unif, model.abm.bc_unif)
    Nbc = _kw(kwargs, :Nbc, model.abm.Nbc)
    bc_probs = _kw(kwargs, :bc_probs, model.abm.bc_probs)
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
    _validate_tmax_vector_constraints(exp.tmax, exp.t_Pass)
    _validate_tmax_length(exp.tmax, n_rep)

    barcode_kwargs = (; skew_lib = skew_lib, use_lib_probs = use_lib_probs,
                       bc_unif = bc_unif, Nbc = Nbc, bc_probs = bc_probs)

    n_pass_eff = exp.tmax isa AbstractVector ? 1 : (length(_passage_times(exp.t_Pass, Float64(exp.tmax))) + 1)

    model_eff = _with_drug_effect(model, de)
    rep_cells = _expand_split_cells_abm(model_eff, exp, n_rep;
                                        R_real = R_real,
                                        drug_effect = de,
                                        skew_lib = skew_lib,
                                        use_lib_probs = use_lib_probs,
                                        split_after_barcoding = split_after_barcoding,
                                        bc_unif = bc_unif,
                                        Nbc = Nbc,
                                        bc_probs = bc_probs,
                                        dt_save_at = dt_save_at,
                                        t_frac = t_frac)

    fin_t_outs = Float64[]
    fin_u_outs = Float64[]
    lin_df_outs = DataFrame[]
    sub_lin_df_outs = DataFrame[]
    sim_dfs = DataFrame[]

    if !just_lin && run_colony
        col_cells_tx1 = seed_cells(nCol, model.params.rho, Int64(1e6); barcode_kwargs...)
        col_cells_tx0 = seed_cells(nCol, model.params.rho, Int64(1e6); barcode_kwargs...)
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
        rep_tmax = _replicate_tmax(exp.tmax, n_rep, i)

        if !just_lin && run_IC
            IC_cells_1 = seed_cells(IC_n0, model.params.rho, model.abm.Nbuff; barcode_kwargs...)
            IC_cells_0 = seed_cells(IC_n0, model.params.rho, model.abm.Nbuff; barcode_kwargs...)

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
            t0 = 0.0, tmax = rep_tmax, t_Pass = exp.t_Pass,
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
    Nvec::Vector{Int64}, nS_vec::Vector{Int64}, nDS_vec::Vector{Int64}, nDR_vec::Vector{Int64}, nR_vec::Vector{Int64},
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

    update_track_vec_resdmg!(kmc_out, Nvec, nS_vec, nDS_vec, nDR_vec, nR_vec, tvec, Pvec)
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
        kmc_out = run_model_core_abm(model, state, sim; treat = treat)
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
    R_real::String = "b",
    drug_effect::Symbol = model.params.drug_effect,
    skew_lib::Bool = model.abm.skew_lib,
    use_lib_probs::Bool = model.abm.use_lib_probs,
    split_after_barcoding::Bool = model.abm.split_after_barcoding,
    bc_unif::Float64 = model.abm.bc_unif,
    Nbc::Int64 = model.abm.Nbc,
    bc_probs::Vector{Float64} = model.abm.bc_probs,
    dt_save_at::Float64 = model.abm.dt_save_at,
    t_frac::Float64 = model.abm.t_frac)

    Nbuff = model.abm.Nbuff
    barcode_kwargs = (; skew_lib = skew_lib, use_lib_probs = use_lib_probs,
                       bc_unif = bc_unif, Nbc = Nbc, bc_probs = bc_probs)

    if split_after_barcoding
        final_seed = _nseed_last(exp.Nseed)
        n_rep * final_seed <= exp.n0 || error("split_after_barcoding requires n0 >= n_rep*Nseed (n0=$(exp.n0), n_rep=$(n_rep), Nseed=$(final_seed)).")

        exp_cells = seed_resdmg_cells(exp.n0, 0.0, Nbuff;
                                      barcode_kwargs...)
        _assign_resistance_by_barcode!(exp_cells, model.params.rho, exp.n0)

        rep_cells = sample(alive_cells(exp_cells), n_rep * final_seed, replace = false)
        rep_cells = reshape(rep_cells, (final_seed, n_rep))

        fin_rep_cells = Vector{Vector{ResDmgCell}}(undef, n_rep)
        for i in 1:n_rep
            fin_rep_cells[i] = collect(rep_cells[:, i])
            extend_with_dead_cells!(fin_rep_cells[i], Nbuff, make_dead_resdmg_cell)
        end
        return fin_rep_cells
    end

    exp_cells = seed_resdmg_cells(exp.n0, model.params.rho, Nbuff;
                                  barcode_kwargs...)

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
    use_lib_probs = _kw(kwargs, :use_lib_probs, model.abm.use_lib_probs)
    split_after_barcoding = _kw(kwargs, :split_after_barcoding, model.abm.split_after_barcoding)
    bc_unif = _kw(kwargs, :bc_unif, model.abm.bc_unif)
    Nbc = _kw(kwargs, :Nbc, model.abm.Nbc)
    bc_probs = _kw(kwargs, :bc_probs, model.abm.bc_probs)
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
    _validate_tmax_vector_constraints(exp.tmax, exp.t_Pass)
    _validate_tmax_length(exp.tmax, n_rep)

    barcode_kwargs = (; skew_lib = skew_lib, use_lib_probs = use_lib_probs,
                       bc_unif = bc_unif, Nbc = Nbc, bc_probs = bc_probs)

    n_pass_eff = exp.tmax isa AbstractVector ? 1 : (length(_passage_times(exp.t_Pass, Float64(exp.tmax))) + 1)

    model_eff = _with_drug_effect(model, de)
    rep_cells = _expand_split_cells_abm(model_eff, exp, n_rep;
                                        R_real = R_real,
                                        drug_effect = de,
                                        skew_lib = skew_lib,
                                        use_lib_probs = use_lib_probs,
                                        split_after_barcoding = split_after_barcoding,
                                        bc_unif = bc_unif,
                                        Nbc = Nbc,
                                        bc_probs = bc_probs,
                                        dt_save_at = dt_save_at,
                                        t_frac = t_frac)

    fin_t_outs = Float64[]
    fin_u_outs = Float64[]
    lin_df_outs = DataFrame[]
    sub_lin_df_outs = DataFrame[]
    sim_dfs = DataFrame[]

    if !just_lin && run_colony
        col_cells_tx1 = seed_resdmg_cells(nCol, model.params.rho, Int64(1e6); barcode_kwargs...)
        col_cells_tx0 = seed_resdmg_cells(nCol, model.params.rho, Int64(1e6); barcode_kwargs...)
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
        rep_tmax = _replicate_tmax(exp.tmax, n_rep, i)

        if !just_lin && run_IC
            IC_cells_1 = seed_resdmg_cells(IC_n0, model.params.rho, model.abm.Nbuff; barcode_kwargs...)
            IC_cells_0 = seed_resdmg_cells(IC_n0, model.params.rho, model.abm.Nbuff; barcode_kwargs...)

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
            t0 = 0.0, tmax = rep_tmax, t_Pass = exp.t_Pass,
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
    K::Int64 = 0
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
    kmc_out = run_model_core_abm(model, state, sim; treat = treat)

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
                         samp_cell_lin_df_vec = samp_cell_lin_df_vec)

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
    K::Int64 = 0
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
    kmc_out = run_model_core_abm(model, state, sim; treat = treat)

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
                                samp_cell_lin_df_vec = samp_cell_lin_df_vec)

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
    R_real = _kw(kwargs, :R_real, "b")
    t_frac = _kw(kwargs, :t_frac, model.abm.t_frac)
    de = normalize_respop_drug_effect(_kw(kwargs, :drug_effect, model.params.drug_effect))
    drug_treatment = _kw(kwargs, :drug_treatment, sim.drug_treatment)
    skew_lib = _kw(kwargs, :skew_lib, model.abm.skew_lib)
    use_lib_probs = _kw(kwargs, :use_lib_probs, model.abm.use_lib_probs)
    bc_unif = _kw(kwargs, :bc_unif, model.abm.bc_unif)
    Nbc = _kw(kwargs, :Nbc, model.abm.Nbc)
    bc_probs = _kw(kwargs, :bc_probs, model.abm.bc_probs)
    dt_save_at = _kw(kwargs, :dt_save_at, model.abm.dt_save_at)

    model_eff = _with_drug_effect(model, de)
    barcode_kwargs = (; skew_lib = skew_lib, use_lib_probs = use_lib_probs,
                       bc_unif = bc_unif, Nbc = Nbc, bc_probs = bc_probs)
    cells = seed_cells(sim.n0, model.params.rho, model.abm.Nbuff;
                       barcode_kwargs...)

    sim_out = _run_abm_simple!(
        model_eff, cells;
        t0 = 0.0, tmax = sim.tmax,
        Nmax = sim.Nmax, Cc = sim.Cc,
        treat_ons = sim.treat_ons, treat_offs = sim.treat_offs,
        dt_save_at = dt_save_at,
        R_real = R_real, t_frac = t_frac, rep = 1,
        treat = drug_treatment, drug_effect = de,
        sub_sample_cells = false, K = 0
    )

    sol_df = DataFrame(
        t = sim_out["tvec"],
        N = sim_out["Nvec"],
        nS = sim_out["nS_vec"],
        nR = sim_out["nR_vec"],
        nE = sim_out["nE_vec"]
    )

    return Dict(
        "lin_df" => join_dfs(sim_out["cell_lin_df_vec"], "bc"),
        "sol_df" => sol_df
    )
end

function _simulate_simple_abm(model::ResDmg_ABM, sim::SimpleSimParams; kwargs...)
    R_real = _kw(kwargs, :R_real, "b")
    t_frac = _kw(kwargs, :t_frac, model.abm.t_frac)
    de = normalize_resdmg_drug_effect(_kw(kwargs, :drug_effect, model.params.drug_effect))
    drug_treatment = _kw(kwargs, :drug_treatment, sim.drug_treatment)
    skew_lib = _kw(kwargs, :skew_lib, model.abm.skew_lib)
    use_lib_probs = _kw(kwargs, :use_lib_probs, model.abm.use_lib_probs)
    bc_unif = _kw(kwargs, :bc_unif, model.abm.bc_unif)
    Nbc = _kw(kwargs, :Nbc, model.abm.Nbc)
    bc_probs = _kw(kwargs, :bc_probs, model.abm.bc_probs)
    dt_save_at = _kw(kwargs, :dt_save_at, model.abm.dt_save_at)

    model_eff = _with_drug_effect(model, de)
    barcode_kwargs = (; skew_lib = skew_lib, use_lib_probs = use_lib_probs,
                       bc_unif = bc_unif, Nbc = Nbc, bc_probs = bc_probs)
    cells = seed_resdmg_cells(sim.n0, model.params.rho, model.abm.Nbuff;
                              barcode_kwargs...)

    sim_out = _run_abm_simple!(
        model_eff, cells;
        t0 = 0.0, tmax = sim.tmax,
        Nmax = sim.Nmax, Cc = sim.Cc,
        treat_ons = sim.treat_ons, treat_offs = sim.treat_offs,
        dt_save_at = dt_save_at,
        R_real = R_real, t_frac = t_frac, rep = 1,
        treat = drug_treatment, drug_effect = de,
        sub_sample_cells = false, K = 0
    )

    sol_df = DataFrame(
        t = sim_out["tvec"],
        N = sim_out["Nvec"],
        nS = sim_out["nS_vec"],
        nDS = sim_out["nDS_vec"],
        nDR = sim_out["nDR_vec"],
        nR = sim_out["nR_vec"]
    )

    return Dict(
        "lin_df" => join_dfs(sim_out["cell_lin_df_vec"], "bc"),
        "sol_df" => sol_df
    )
end

function _expand_split_cells_abm(model::ResPopInVivo_ABM, exp::ExperimentParams, n_rep::Int64;
    R_real::String = "b",
    drug_effect::Symbol = model.params.drug_effect,
    skew_lib::Bool = model.abm.skew_lib,
    use_lib_probs::Bool = model.abm.use_lib_probs,
    split_after_barcoding::Bool = model.abm.split_after_barcoding,
    bc_unif::Float64 = model.abm.bc_unif,
    Nbc::Int64 = model.abm.Nbc,
    bc_probs::Vector{Float64} = model.abm.bc_probs,
    dt_save_at::Float64 = model.abm.dt_save_at,
    t_frac::Float64 = model.abm.t_frac,
    rep_design = _experiment_condition_design(true, false, n_rep),
    inc_pot::Bool = false,
    pot_outputs = nothing)

    Nbuff = model.abm.Nbuff
    barcode_kwargs = (; skew_lib = skew_lib, use_lib_probs = use_lib_probs,
                       bc_unif = bc_unif, Nbc = Nbc, bc_probs = bc_probs)

    if split_after_barcoding
        final_seed = _nseed_last(exp.Nseed)
        n_batches = length(rep_design)
        n_batches * final_seed <= exp.n0 || error("split_after_barcoding requires n0 >= n_batches*Nseed (n0=$(exp.n0), n_batches=$(n_batches), Nseed=$(final_seed)).")

        exp_cells = seed_invivo_cells(exp.n0, 0.0, 0.0, Nbuff;
                                      barcode_kwargs...)
        _assign_resistance_by_barcode!(exp_cells, model.params.rho, exp.n0)
        _assign_engraftment_by_barcode!(exp_cells, model.params.fEG1, exp.n0)

        if inc_pot && pot_outputs !== nothing
            pot_sim = Dict(
                "Nvec" => Int64[length(alive_cells(exp_cells))],
                "tvec" => Float64[0.0],
                "Pvec" => Int64[0],
                "nS_vec" => Int64[sum(cell -> cell.alive && !cell.R && !cell.E, exp_cells)],
                "nR_vec" => Int64[sum(cell -> cell.alive && cell.R && !cell.E, exp_cells)],
                "nE_vec" => Int64[sum(cell -> cell.alive && cell.E, exp_cells)],
                "nS_EG0_vec" => Int64[sum(cell -> cell.alive && !cell.R && !cell.E && !cell.EG, exp_cells)],
                "nS_EG1_vec" => Int64[sum(cell -> cell.alive && !cell.R && !cell.E && cell.EG, exp_cells)],
                "nR_EG0_vec" => Int64[sum(cell -> cell.alive && cell.R && !cell.E && !cell.EG, exp_cells)],
                "nR_EG1_vec" => Int64[sum(cell -> cell.alive && cell.R && !cell.E && cell.EG, exp_cells)],
                "nE_EG0_vec" => Int64[sum(cell -> cell.alive && cell.E && !cell.EG, exp_cells)],
                "nE_EG1_vec" => Int64[sum(cell -> cell.alive && cell.E && cell.EG, exp_cells)]
            )
            pot_outputs["sol_df"] = _invivo_abm_sol_df(pot_sim; cond = "POT", rep = 0)
            pot_outputs["lin_df"] = get_counts(alive_cells(exp_cells), "POT_P0")
        end

        rep_cells = sample(alive_cells(exp_cells), n_batches * final_seed, replace = false)
        rep_cells = reshape(rep_cells, (final_seed, n_batches))

        fin_rep_cells = Vector{Vector{InVivoCancerCell}}(undef, n_batches)
        engraft_rows = DataFrame[]
        for i in 1:n_batches
            design = rep_design[i]
            cells_i = collect(rep_cells[:, i])
            cells_i, stats = engraftment_selection(cells_i, model.params.pEG, model.params.sEG)
            push!(engraft_rows, DataFrame(cond = design.cond,
                                          rep = design.rep,
                                          passage = 1,
                                          N_engraft = stats["N_engraft"],
                                          nEG0_engraft = stats["nEG0_engraft"],
                                          nEG1_engraft = stats["nEG1_engraft"]))
            extend_with_dead_cells!(cells_i, Nbuff, make_dead_cell_invivo)
            fin_rep_cells[i] = cells_i
        end

        engraft_df = isempty(engraft_rows) ? DataFrame(cond = String[], rep = Int[], passage = Int[], N_engraft = Int[], nEG0_engraft = Int[], nEG1_engraft = Int[]) : vcat(engraft_rows...)
        return fin_rep_cells, engraft_df
    end

    exp_cells = seed_invivo_cells(exp.n0, model.params.rho, model.params.fEG1, Nbuff;
                                  barcode_kwargs...)

    expansion_model = ResPopInVivo_ABM(_copy_respop_invivo_params(model.params; al = 0.0, drug_effect = drug_effect);
                                       abm = model.abm)
    pot_sim = Dict(
        "Nvec" => Int64[],
        "tvec" => Float64[],
        "Pvec" => Int64[],
        "nS_vec" => Int64[],
        "nR_vec" => Int64[],
        "nE_vec" => Int64[],
        "nS_EG0_vec" => Int64[],
        "nS_EG1_vec" => Int64[],
        "nR_EG0_vec" => Int64[],
        "nR_EG1_vec" => Int64[],
        "nE_EG0_vec" => Int64[],
        "nE_EG1_vec" => Int64[]
    )
    pot_t_offset = 0.0

    if (exp.t_exp isa AbstractVector) && (exp.Nseed isa AbstractVector)
        @assert length(exp.t_exp) == length(exp.Nseed) "t_exp and Nseed vectors must be of same length"

        for i in 1:(length(exp.t_exp) - 1)
            stage_sim = ABMSimParams(
                t0 = 0.0,
                tmax = Float64(exp.t_exp[i]),
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
            kmc_out = run_model_core_abm(expansion_model, ResPopInVivoABMState(exp_cells), stage_sim; treat = false)
            if inc_pot
                _append_invivo_abm_solution_outputs!(kmc_out,
                                                     pot_sim["Nvec"], pot_sim["nS_vec"], pot_sim["nR_vec"], pot_sim["nE_vec"],
                                                     pot_sim["tvec"], pot_sim["Pvec"],
                                                     pot_sim["nS_EG0_vec"], pot_sim["nS_EG1_vec"],
                                                     pot_sim["nR_EG0_vec"], pot_sim["nR_EG1_vec"],
                                                     pot_sim["nE_EG0_vec"], pot_sim["nE_EG1_vec"];
                                                     t_offset = pot_t_offset, passage = 0)
                pot_t_offset += last(kmc_out.tvec)
            end

            exp_cells = alive_cells(exp_cells)
            length(exp_cells) >= exp.Nseed[i] || error("Not enough cells after expansion at stage $i for bottlenecking.")
            exp_cells = sample(exp_cells, Int64(exp.Nseed[i]), replace = false)
            extend_with_dead_cells!(exp_cells, Nbuff, make_dead_cell_invivo)
        end

        final_sim = ABMSimParams(
            t0 = 0.0,
            tmax = Float64(exp.t_exp[end]),
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
        kmc_out = run_model_core_abm(expansion_model, ResPopInVivoABMState(exp_cells), final_sim; treat = false)
        if inc_pot
            _append_invivo_abm_solution_outputs!(kmc_out,
                                                 pot_sim["Nvec"], pot_sim["nS_vec"], pot_sim["nR_vec"], pot_sim["nE_vec"],
                                                 pot_sim["tvec"], pot_sim["Pvec"],
                                                 pot_sim["nS_EG0_vec"], pot_sim["nS_EG1_vec"],
                                                 pot_sim["nR_EG0_vec"], pot_sim["nR_EG1_vec"],
                                                 pot_sim["nE_EG0_vec"], pot_sim["nE_EG1_vec"];
                                                 t_offset = pot_t_offset, passage = 0)
        end
        exp_cells = alive_cells(exp_cells)
        final_seed = Int64(exp.Nseed[end])
    elseif (exp.t_exp isa Real) && (exp.Nseed isa Integer)
        final_sim = ABMSimParams(
            t0 = 0.0,
            tmax = Float64(exp.t_exp),
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
        kmc_out = run_model_core_abm(expansion_model, ResPopInVivoABMState(exp_cells), final_sim; treat = false)
        if inc_pot
            _append_invivo_abm_solution_outputs!(kmc_out,
                                                 pot_sim["Nvec"], pot_sim["nS_vec"], pot_sim["nR_vec"], pot_sim["nE_vec"],
                                                 pot_sim["tvec"], pot_sim["Pvec"],
                                                 pot_sim["nS_EG0_vec"], pot_sim["nS_EG1_vec"],
                                                 pot_sim["nR_EG0_vec"], pot_sim["nR_EG1_vec"],
                                                 pot_sim["nE_EG0_vec"], pot_sim["nE_EG1_vec"];
                                                 t_offset = pot_t_offset, passage = 0)
        end
        exp_cells = alive_cells(exp_cells)
        final_seed = Int64(exp.Nseed)
    else
        error("t_exp and Nseed must both be scalars or both be vectors of equal length.")
    end

    if inc_pot && pot_outputs !== nothing
        pot_outputs["sol_df"] = _invivo_abm_sol_df(pot_sim; cond = "POT", rep = 0)
        pot_outputs["lin_df"] = get_counts(exp_cells, "POT_P0")
    end

    n_batches = length(rep_design)
    n_batches * final_seed <= length(exp_cells) || error("Not enough cells for $n_batches replicates of size $final_seed.")
    rep_cells = sample(exp_cells, n_batches * final_seed, replace = false)
    rep_cells = reshape(rep_cells, (final_seed, n_batches))

    fin_rep_cells = Vector{Vector{InVivoCancerCell}}(undef, n_batches)
    engraft_rows = DataFrame[]
    for i in 1:n_batches
        design = rep_design[i]
        cells_i = collect(rep_cells[:, i])
        cells_i, stats = engraftment_selection(cells_i, model.params.pEG, model.params.sEG)
        push!(engraft_rows, DataFrame(cond = design.cond,
                                      rep = design.rep,
                                      passage = 1,
                                      N_engraft = stats["N_engraft"],
                                      nEG0_engraft = stats["nEG0_engraft"],
                                      nEG1_engraft = stats["nEG1_engraft"]))
        extend_with_dead_cells!(cells_i, Nbuff, make_dead_cell_invivo)
        fin_rep_cells[i] = cells_i
    end

    engraft_df = isempty(engraft_rows) ? DataFrame(cond = String[], rep = Int[], passage = Int[], N_engraft = Int[], nEG0_engraft = Int[], nEG1_engraft = Int[]) : vcat(engraft_rows...)
    return fin_rep_cells, engraft_df
end

function _run_abm_passage_experiment_invivo!(
    model::ResPopInVivo_ABM,
    cells::Vector{InVivoCancerCell};
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
    cond::String = "DT",
    sub_sample_cells::Bool = false,
    K::Int64 = 0
)
    t_pass_vec = _passage_times(t_Pass, tmax)
    boundaries = vcat([t0], filter(x -> x < tmax, t_pass_vec), [tmax])

    cell_lin_df_vec = DataFrame[]
    samp_cell_lin_df_vec = DataFrame[]
    Nvec = Int64[]
    nS_vec = Int64[]
    nR_vec = Int64[]
    nE_vec = Int64[]
    nS_EG0_vec = Int64[]
    nS_EG1_vec = Int64[]
    nR_EG0_vec = Int64[]
    nR_EG1_vec = Int64[]
    nE_EG0_vec = Int64[]
    nE_EG1_vec = Int64[]
    tvec = Float64[]
    Pvec = Int64[]
    engraft_rows = DataFrame[]

    for seg_idx in 1:(length(boundaries) - 1)
        seg_sim = ABMSimParams(
            t0 = boundaries[seg_idx],
            tmax = boundaries[seg_idx + 1],
            Nmax = Nmax,
            Cc = Cc,
            treat_ons = treat_ons,
            treat_offs = treat_offs,
            dt_save_at = dt_save_at,
            R_real = R_real,
            t_frac = t_frac,
            Passage = seg_idx,
            drug_effect = drug_effect
        )
        kmc_out = run_model_core_abm(model, ResPopInVivoABMState(cells), seg_sim; treat = treat)

        _record_abm_outputs!(kmc_out, cells, rep, seg_idx, cell_lin_df_vec,
                             Nvec, nS_vec, nR_vec, nE_vec, tvec, Pvec,
                             sub_sample_cells = sub_sample_cells, K = K,
                             samp_cell_lin_df_vec = samp_cell_lin_df_vec,
                             cond = cond)
        _append_invivo_abm_eg_outputs!(kmc_out,
                                       nS_EG0_vec, nS_EG1_vec,
                                       nR_EG0_vec, nR_EG1_vec,
                                       nE_EG0_vec, nE_EG1_vec)

        if seg_idx < (length(boundaries) - 1)
            live_cells = alive_cells(cells)
            if length(live_cells) < Nseed
                break
            end
            live_cells = sample(live_cells, Nseed, replace = false)
            live_cells, stats = engraftment_selection(live_cells, model.params.pEG, model.params.sEG)
            push!(engraft_rows, DataFrame(cond = cond,
                                          rep = rep,
                                          passage = seg_idx + 1,
                                          N_engraft = stats["N_engraft"],
                                          nEG0_engraft = stats["nEG0_engraft"],
                                          nEG1_engraft = stats["nEG1_engraft"]))
            cells = live_cells
            extend_with_dead_cells!(cells, Nbuff, make_dead_cell_invivo)
        end
    end

    out = Dict(
        "cell_lin_df_vec" => cell_lin_df_vec,
        "Nvec" => Nvec,
        "tvec" => tvec,
        "Pvec" => Pvec,
        "nS_vec" => nS_vec,
        "nR_vec" => nR_vec,
        "nE_vec" => nE_vec,
        "nS_EG0_vec" => nS_EG0_vec,
        "nS_EG1_vec" => nS_EG1_vec,
        "nR_EG0_vec" => nR_EG0_vec,
        "nR_EG1_vec" => nR_EG1_vec,
        "nE_EG0_vec" => nE_EG0_vec,
        "nE_EG1_vec" => nE_EG1_vec,
        "engraft_df" => isempty(engraft_rows) ? DataFrame(cond = String[], rep = Int[], passage = Int[], N_engraft = Int[], nEG0_engraft = Int[], nEG1_engraft = Int[]) : vcat(engraft_rows...)
    )
    if sub_sample_cells
        out["sub_samp_cell_lin_df_vec"] = samp_cell_lin_df_vec
    end
    return out
end

function _simulate_experiment_abm(model::ResPopInVivo_ABM, exp::ExperimentParams; kwargs...)
    n_rep = _kw(kwargs, :n_rep, exp.n_rep)
    R_real = _kw(kwargs, :R_real, "b")
    t_frac = _kw(kwargs, :t_frac, model.abm.t_frac)
    just_lin = _kw(kwargs, :just_lin, false)
    de = normalize_respop_drug_effect(_kw(kwargs, :drug_effect, model.params.drug_effect))
    drug_treatment = _kw(kwargs, :drug_treatment, exp.drug_treatment)
    inc_control = _kw(kwargs, :inc_control, exp.inc_control)
    inc_pot = _kw(kwargs, :inc_pot, exp.inc_pot)
    sub_sample_cells = _kw(kwargs, :sub_sample_cells, model.abm.sub_sample_cells)
    K = _kw(kwargs, :K, model.abm.K)
    skew_lib = _kw(kwargs, :skew_lib, model.abm.skew_lib)
    use_lib_probs = _kw(kwargs, :use_lib_probs, model.abm.use_lib_probs)
    split_after_barcoding = _kw(kwargs, :split_after_barcoding, model.abm.split_after_barcoding)
    bc_unif = _kw(kwargs, :bc_unif, model.abm.bc_unif)
    Nbc = _kw(kwargs, :Nbc, model.abm.Nbc)
    bc_probs = _kw(kwargs, :bc_probs, model.abm.bc_probs)
    dt_save_at = _kw(kwargs, :dt_save_at, model.abm.dt_save_at)

    _validate_tmax_vector_constraints(exp.tmax, exp.t_Pass)
    _validate_tmax_length(exp.tmax, n_rep)

    model_eff = _with_drug_effect(model, de)
    rep_design = _experiment_condition_design(drug_treatment, inc_control, n_rep)
    pot_outputs = inc_pot ? Dict{String, Any}() : nothing
    rep_cells, split_engraft_df = _expand_split_cells_abm(model_eff, exp, n_rep;
                                                          R_real = R_real,
                                                          drug_effect = de,
                                                          skew_lib = skew_lib,
                                                          use_lib_probs = use_lib_probs,
                                                          split_after_barcoding = split_after_barcoding,
                                                          bc_unif = bc_unif,
                                                          Nbc = Nbc,
                                                          bc_probs = bc_probs,
                                                          dt_save_at = dt_save_at,
                                                          t_frac = t_frac,
                                                          rep_design = rep_design,
                                                          inc_pot = inc_pot,
                                                          pot_outputs = pot_outputs)

    fin_t_outs = Float64[]
    fin_u_outs = Float64[]
    fin_cond_outs = String[]
    fin_rep_outs = Int64[]
    lin_df_outs = DataFrame[]
    sub_lin_df_outs = DataFrame[]
    sim_dfs = DataFrame[]
    engraft_rows = DataFrame[split_engraft_df]

    if inc_pot
        push!(lin_df_outs, pot_outputs["lin_df"])
        if !just_lin
            push!(sim_dfs, pot_outputs["sol_df"])
        end
    end

    nseed_last = _nseed_last(exp.Nseed)
    for i in eachindex(rep_design)
        design = rep_design[i]
        rep_tmax = _replicate_tmax(exp.tmax, n_rep, design.rep)

        extend_with_dead_cells!(rep_cells[i], model.abm.Nbuff, make_dead_cell_invivo)
        sim = _run_abm_passage_experiment_invivo!(
            model_eff, rep_cells[i];
            t0 = 0.0, tmax = rep_tmax, t_Pass = exp.t_Pass,
            Nseed = nseed_last, Nmax = exp.Nmax, Cc = exp.Cc,
            treat_ons = exp.treat_ons, treat_offs = exp.treat_offs,
            dt_save_at = dt_save_at, Nbuff = model.abm.Nbuff,
            R_real = R_real, t_frac = t_frac, rep = design.rep,
            treat = design.treat, drug_effect = de,
            cond = design.cond,
            sub_sample_cells = sub_sample_cells, K = K
        )

        push!(engraft_rows, sim["engraft_df"])
        push!(lin_df_outs, join_dfs(sim["cell_lin_df_vec"], "bc"))
        if sub_sample_cells
            push!(sub_lin_df_outs, join_dfs(sim["sub_samp_cell_lin_df_vec"], "bc"))
        end

        if !just_lin
            push!(sim_dfs, _invivo_abm_sol_df(sim; cond = design.cond, rep = design.rep))
            if !isempty(sim["tvec"])
                push!(fin_t_outs, last(sim["tvec"]))
                push!(fin_u_outs, last(sim["Nvec"]))
                push!(fin_cond_outs, design.cond)
                push!(fin_rep_outs, design.rep)
            end
        end
    end

    out = Dict{String, Any}(
        "lin_df" => join_dfs(lin_df_outs, "bc"),
        "engraft_df" => isempty(engraft_rows) ? DataFrame(cond = String[], rep = Int[], passage = Int[], N_engraft = Int[], nEG0_engraft = Int[], nEG1_engraft = Int[]) : vcat(engraft_rows...)
    )

    if sub_sample_cells
        out["sub_lin_df"] = join_dfs(sub_lin_df_outs, "bc")
    end

    if !just_lin
        out["t"] = fin_t_outs
        out["u"] = fin_u_outs
        out["cond"] = fin_cond_outs
        out["rep"] = fin_rep_outs
        out["sol_df"] = isempty(sim_dfs) ? DataFrame() : vcat(sim_dfs...)
    end

    return out
end

function _run_abm_simple!(
    model::ResPopInVivo_ABM,
    cells::Vector{InVivoCancerCell};
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
    K::Int64 = 0
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
    kmc_out = run_model_core_abm(model, ResPopInVivoABMState(cells), sim; treat = treat)

    cell_lin_df_vec = DataFrame[]
    samp_cell_lin_df_vec = DataFrame[]
    Nvec = Int64[]
    nS_vec = Int64[]
    nR_vec = Int64[]
    nE_vec = Int64[]
    nS_EG0_vec = Int64[]
    nS_EG1_vec = Int64[]
    nR_EG0_vec = Int64[]
    nR_EG1_vec = Int64[]
    nE_EG0_vec = Int64[]
    nE_EG1_vec = Int64[]
    tvec = Float64[]
    Pvec = Int64[]

    _record_abm_outputs!(kmc_out, cells, rep, 1, cell_lin_df_vec,
                         Nvec, nS_vec, nR_vec, nE_vec, tvec, Pvec,
                         sub_sample_cells = sub_sample_cells, K = K,
                         samp_cell_lin_df_vec = samp_cell_lin_df_vec)
    _append_invivo_abm_eg_outputs!(kmc_out,
                                   nS_EG0_vec, nS_EG1_vec,
                                   nR_EG0_vec, nR_EG1_vec,
                                   nE_EG0_vec, nE_EG1_vec)

    out = Dict(
        "cell_lin_df_vec" => cell_lin_df_vec,
        "Nvec" => Nvec,
        "tvec" => tvec,
        "Pvec" => Pvec,
        "nS_vec" => nS_vec,
        "nR_vec" => nR_vec,
        "nE_vec" => nE_vec,
        "nS_EG0_vec" => nS_EG0_vec,
        "nS_EG1_vec" => nS_EG1_vec,
        "nR_EG0_vec" => nR_EG0_vec,
        "nR_EG1_vec" => nR_EG1_vec,
        "nE_EG0_vec" => nE_EG0_vec,
        "nE_EG1_vec" => nE_EG1_vec
    )
    if sub_sample_cells
        out["sub_samp_cell_lin_df_vec"] = samp_cell_lin_df_vec
    end
    return out
end

function _simulate_simple_abm(model::ResPopInVivo_ABM, sim::SimpleSimParams; kwargs...)
    R_real = _kw(kwargs, :R_real, "b")
    t_frac = _kw(kwargs, :t_frac, model.abm.t_frac)
    de = normalize_respop_drug_effect(_kw(kwargs, :drug_effect, model.params.drug_effect))
    drug_treatment = _kw(kwargs, :drug_treatment, sim.drug_treatment)
    skew_lib = _kw(kwargs, :skew_lib, model.abm.skew_lib)
    use_lib_probs = _kw(kwargs, :use_lib_probs, model.abm.use_lib_probs)
    bc_unif = _kw(kwargs, :bc_unif, model.abm.bc_unif)
    Nbc = _kw(kwargs, :Nbc, model.abm.Nbc)
    bc_probs = _kw(kwargs, :bc_probs, model.abm.bc_probs)
    dt_save_at = _kw(kwargs, :dt_save_at, model.abm.dt_save_at)

    model_eff = _with_drug_effect(model, de)
    barcode_kwargs = (; skew_lib = skew_lib, use_lib_probs = use_lib_probs,
                       bc_unif = bc_unif, Nbc = Nbc, bc_probs = bc_probs)
    cells = seed_invivo_cells(sim.n0, model.params.rho, model.params.fEG1, model.abm.Nbuff;
                              barcode_kwargs...)

    sim_out = _run_abm_simple!(
        model_eff, cells;
        t0 = 0.0, tmax = sim.tmax,
        Nmax = sim.Nmax, Cc = sim.Cc,
        treat_ons = sim.treat_ons, treat_offs = sim.treat_offs,
        dt_save_at = dt_save_at,
        R_real = R_real, t_frac = t_frac, rep = 1,
        treat = drug_treatment, drug_effect = de,
        sub_sample_cells = false, K = 0
    )

    sol_df = DataFrame(
        t = sim_out["tvec"],
        nS_EG0 = sim_out["nS_EG0_vec"],
        nS_EG1 = sim_out["nS_EG1_vec"],
        nR_EG0 = sim_out["nR_EG0_vec"],
        nR_EG1 = sim_out["nR_EG1_vec"],
        nE_EG0 = sim_out["nE_EG0_vec"],
        nE_EG1 = sim_out["nE_EG1_vec"],
        nS = sim_out["nS_vec"],
        nR = sim_out["nR_vec"],
        nE = sim_out["nE_vec"],
        n_EG0 = sim_out["nS_EG0_vec"] .+ sim_out["nR_EG0_vec"] .+ sim_out["nE_EG0_vec"],
        n_EG1 = sim_out["nS_EG1_vec"] .+ sim_out["nR_EG1_vec"] .+ sim_out["nE_EG1_vec"],
        N = sim_out["Nvec"]
    )

    return Dict(
        "lin_df" => join_dfs(sim_out["cell_lin_df_vec"], "bc"),
        "sol_df" => sol_df
    )
end

