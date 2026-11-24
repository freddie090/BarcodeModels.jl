"""Construct the optional fixed-interval, top-barcode trajectory recorder."""
function _lin_track_recorder(model, settings, kwargs)
    enabled = _kw(kwargs, :lin_track, settings.lin_track)
    rec_t = _kw(kwargs, :lin_track_rec_t, settings.lin_track_rec_t)
    top_x = _kw(kwargs, :top_x, settings.top_x)
    validate_lin_track_params(enabled, rec_t, top_x)
    enabled || return nothing
    df = DataFrame(snapshot = Int64[], t = Float64[], cond = String[], rep = Int64[],
        passage = Int64[], event = String[], bc = Union{Missing, Float64}[],
        N = Int64[], population_N = Int64[], nS = Int64[], nR = Int64[])
    if model isa Union{ResDmg_ABM, ResDmg_ABM_EvBC}
        df.nDS = Int64[]
        df.nDR = Int64[]
    else
        df.nE = Int64[]
    end
    return LinTrackRecorder(Float64(rec_t), Int64(top_x), df, "DT", 1, 0.0, 1, 0)
end

_lin_track_context!(::Nothing, cond, rep) = nothing
function _lin_track_context!(tracker::LinTrackRecorder, cond, rep)
    tracker.cond = cond
    tracker.rep = rep
    tracker.offset = 0.0
    tracker.next_tick = 1
    return nothing
end

_lin_track_pheno(cell::Union{CancerCell, CancerCellEvBC, InVivoCancerCell}) =
    cell.E ? 3 : cell.R ? 2 : 1
_lin_track_pheno(cell::Union{ResDmgCell, ResDmgCellEvBC}) =
    cell.DS ? 3 : cell.DR ? 4 : cell.R ? 2 : 1

_lin_track_snapshot!(::Nothing, cells, time, passage, event) = nothing
function _lin_track_snapshot!(tracker::LinTrackRecorder, cells, time, passage, event)
    counts = Dict{Float64, Vector{Int64}}()
    population_N = Int64(0)
    for cell in cells
        cell.alive || continue
        population_N += 1
        abundance = get!(counts, cell.barcode) do
            zeros(Int64, 4)
        end
        abundance[_lin_track_pheno(cell)] += 1
    end
    # Stable barcode tie-breaking consumes no simulation random numbers.
    barcodes = sort!(collect(keys(counts)); by = bc -> (-sum(counts[bc]), bc))
    resize!(barcodes, min(length(barcodes), tracker.top_x))
    tracker.snapshot += 1
    # A missing barcode marks an empty population, including complete extinction.
    for bc in (isempty(barcodes) ? [missing] : barcodes)
        n = ismissing(bc) ? zeros(Int64, 4) : counts[bc]
        row = (snapshot = tracker.snapshot, t = Float64(time), cond = tracker.cond,
            rep = tracker.rep, passage = tracker.cond == "POT" ? 0 : passage,
            event = event, bc = bc, N = sum(n), population_N = population_N, nS = n[1], nR = n[2])
        phenotypes = hasproperty(tracker.counts, :nE) ? (nE = n[3],) : (nDS = n[3], nDR = n[4])
        push!(tracker.counts, merge(row, phenotypes))
    end
    return nothing
end

_lin_track_until!(::Nothing, cells, time, passage) = nothing
function _lin_track_until!(tracker::LinTrackRecorder, cells, time, passage)
    # Evaluate the unchanged pre-event state at every crossed observation time.
    stop = tracker.offset + time
    while tracker.next_tick * tracker.rec_t <= stop ||
          isapprox(tracker.next_tick * tracker.rec_t, stop; rtol = 8eps(Float64), atol = 0.0)
        tick = min(tracker.next_tick * tracker.rec_t, stop)
        _lin_track_snapshot!(tracker, cells, tick, passage, "grid")
        tracker.next_tick += 1
    end
    return nothing
end

function _lin_track_start!(sim::ABMSimParams, cells)
    tracker = sim.lin_tracker
    tracker === nothing && return nothing
    _lin_track_snapshot!(tracker, cells, tracker.offset + sim.t0, sim.Passage, "start")
    return nothing
end

function _lin_track_end!(sim::ABMSimParams, cells, time, out)
    tracker = sim.lin_tracker
    tracker === nothing && return out
    time = min(time, sim.tmax)
    _lin_track_until!(tracker, cells, time, sim.Passage)
    _lin_track_snapshot!(tracker, cells, tracker.offset + time, sim.Passage, "end")
    # Expansion stages use local simulation clocks; retain one continuous POT clock.
    tracker.cond == "POT" && (tracker.offset += time - sim.t0)
    return out
end

"""Replace missing entries in an ABM output table with zeros."""
function rep_missings(df::DataFrame)
    for col in names(df)
        df[ismissing.(df[:, col]), col] .= 0.0
    end
    return df
end

"""Join ABM barcode count tables on the requested column."""
function join_dfs(dfs::Vector{DataFrame}, colname::String)
    if length(dfs) == 1
        full_df = dfs[1]
    else
        full_df = dfs[1]
        for i in 2:length(dfs)
            full_df = outerjoin(full_df, dfs[i], on = Symbol(colname))
        end
        full_df = rep_missings(full_df)
    end
    return full_df
end

"""Count barcode abundances in the supplied cells."""
function get_counts(cells, rep_name::String)
    if length(cells) == 0
        df = DataFrame(bc = Float64[], n = Int64[])
    else
        df = DataFrame(bc = map(x -> x.barcode, cells))
        df = combine(groupby(df, :bc), nrow)
    end

    rename!(df, [:bc, Symbol(rep_name)])
    return df
end

"""
    neutral_mutation_counts(cells, population_name="population")

Count heritable neutral mutation IDs among live cells, returning `mutation_id`,
`<population_name>_count` and `<population_name>_frequency` columns. Frequency is
the fraction of live cells carrying the mutation (not allele frequency). Dead
cells are ignored; empty or mutation-free populations return a typed empty table.
Count only supplied cells to summarise a sample. Mutations absent from it are omitted.
"""
function neutral_mutation_counts(cells, population_name::AbstractString = "population")
    counts = Dict{Int64, Int64}()
    population_size = 0
    for cell in cells
        cell.alive || continue
        population_size += 1
        muts = cell.muts
        muts === nothing && continue
        # Simulator mutation vectors contain each inherited infinite-sites ID once.
        for id in muts
            counts[id] = get(counts, id, Int64(0)) + 1
        end
    end
    ids = sort!(collect(keys(counts)))
    abundance = Int64[counts[id] for id in ids]
    frequency = population_size == 0 ? Float64[] : abundance ./ population_size
    return DataFrame(:mutation_id => ids,
        Symbol(population_name, "_count") => abundance,
        Symbol(population_name, "_frequency") => frequency)
end

function _add_neutral_outputs!(out, tracker::NeutralMutationTracker)
    out["mut_df"] = isempty(tracker.counts) ? DataFrame(mutation_id = Int64[]) :
        sort!(join_dfs(tracker.counts, "mutation_id"), :mutation_id)
    if !isempty(tracker.sample_counts)
        out["sub_mut_df"] = sort!(join_dfs(tracker.sample_counts, "mutation_id"), :mutation_id)
    end
    return out
end

_add_neutral_outputs!(out, ::Nothing) = out

_collect_pre_treatment_mutations!(::Nothing, event, rep) = nothing

function _collect_pre_treatment_mutations!(tables::Vector{DataFrame}, event::InstantKillState, rep)
    event.pre_treatment_mutations === nothing && return nothing
    push!(tables, rename(event.pre_treatment_mutations,
        :population_count => Symbol("DT", rep, "_count"),
        :population_frequency => Symbol("DT", rep, "_frequency")))
    return nothing
end

_add_pre_treatment_mutation_outputs!(out, ::Nothing) = out

function _add_pre_treatment_mutation_outputs!(out, tables::Vector{DataFrame})
    out["mut_df_pre_treatment"] = isempty(tables) ? DataFrame(mutation_id = Int64[]) :
        sort!(join_dfs(tables, "mutation_id"), :mutation_id)
    return out
end

function _add_pre_treatment_mutation_outputs!(out, event::InstantKillState)
    event.settings.save_pre_treatment_mutations || return out
    tables = DataFrame[]
    _collect_pre_treatment_mutations!(tables, event, 1)
    return _add_pre_treatment_mutation_outputs!(out, tables)
end

function _instant_kill_row(event::InstantKillState, rep::Int64)
    population_trigger = _instant_population_trigger(event.settings)
    time = population_trigger ? event.time : event.settings.inst_kill_time
    row = (rep = rep, time = time,
        mode = event.settings.inst_kill_mode,
        status = event.applied ? "applied" : "not_reached",
        mutation_id = event.mutation_id, carrier_fraction = event.carrier_fraction,
        N_before = event.N_before, N_killed = event.N_killed, N_after = event.N_after)
    return population_trigger ? merge(row, (inst_kill_N = event.settings.inst_kill_N,)) : row
end

"""Append ResPop population and passage trajectories to output vectors."""
function update_track_vec!(kmc_out,
    Nvec::Vector{Int64},
    nS_vec::Vector{Int64}, nR_vec::Vector{Int64}, nE_vec::Vector{Int64},
    tvec::Vector{Float64}, Pvec::Vector{Int64})

    append!(Nvec, kmc_out.Nvec)
    append!(nS_vec, kmc_out.Svec)
    append!(nR_vec, kmc_out.Rvec)
    append!(nE_vec, kmc_out.Evec)
    append!(tvec, kmc_out.tvec)
    append!(Pvec, kmc_out.Pvec)
end

"""Append ResDmg population and passage trajectories to output vectors."""
function update_track_vec_resdmg!(kmc_out,
    Nvec::Vector{Int64},
    nS_vec::Vector{Int64}, nDS_vec::Vector{Int64}, nDR_vec::Vector{Int64}, nR_vec::Vector{Int64},
    tvec::Vector{Float64}, Pvec::Vector{Int64})

    append!(Nvec, kmc_out.Nvec)
    append!(nS_vec, kmc_out.Svec)
    append!(nDS_vec, kmc_out.DSvec)
    append!(nDR_vec, kmc_out.DRvec)
    append!(nR_vec, kmc_out.Rvec)
    append!(tvec, kmc_out.tvec)
    append!(Pvec, kmc_out.Pvec)
end

"""Record ResPop ABM trajectory, barcode and optional subsample outputs."""
function _record_abm_outputs!(kmc_out, cells, rep::Int64, curr_P::Int64,
    cell_lin_df_vec::Vector{DataFrame},
    Nvec::Vector{Int64}, nS_vec::Vector{Int64}, nR_vec::Vector{Int64}, nE_vec::Vector{Int64},
    tvec::Vector{Float64}, Pvec::Vector{Int64};
    sub_sample_cells::Bool = false, K::Int64 = 0,
    samp_cell_lin_df_vec::Vector{DataFrame} = DataFrame[],
    neutral_tracker::Union{Nothing, NeutralMutationTracker} = nothing)

    live_cells = alive_cells(cells)
    bc_df = get_counts(live_cells, string("DT", rep, "_P", curr_P))
    push!(cell_lin_df_vec, bc_df)

    if neutral_tracker !== nothing
        push!(neutral_tracker.counts, neutral_mutation_counts(live_cells, string("DT", rep, "_P", curr_P)))
    end

    if sub_sample_cells
        if K <= length(live_cells)
            samp_cells = sample(live_cells, K, replace = false)
            push!(samp_cell_lin_df_vec, get_counts(samp_cells, string("DT", rep, "_P", curr_P)))
        else
            push!(samp_cell_lin_df_vec, get_counts(live_cells, string("DT", rep, "_P", curr_P)))
        end
        if neutral_tracker !== nothing
            sampled = K <= length(live_cells) ? samp_cells : live_cells
            push!(neutral_tracker.sample_counts, neutral_mutation_counts(sampled, string("DT", rep, "_P", curr_P)))
        end
    end

    update_track_vec!(kmc_out, Nvec, nS_vec, nR_vec, nE_vec, tvec, Pvec)
end

"""Record ResDmg ABM trajectory, barcode and optional subsample outputs."""
function _record_resdmg_abm_outputs!(kmc_out, cells, rep::Int64, curr_P::Int64,
    cell_lin_df_vec::Vector{DataFrame},
    Nvec::Vector{Int64}, nS_vec::Vector{Int64}, nDS_vec::Vector{Int64}, nDR_vec::Vector{Int64}, nR_vec::Vector{Int64},
    tvec::Vector{Float64}, Pvec::Vector{Int64};
    sub_sample_cells::Bool = false, K::Int64 = 0,
    samp_cell_lin_df_vec::Vector{DataFrame} = DataFrame[],
    neutral_tracker::Union{Nothing, NeutralMutationTracker} = nothing)

    live_cells = alive_cells(cells)
    bc_df = get_counts(live_cells, string("DT", rep, "_P", curr_P))
    push!(cell_lin_df_vec, bc_df)

    if neutral_tracker !== nothing
        push!(neutral_tracker.counts, neutral_mutation_counts(live_cells, string("DT", rep, "_P", curr_P)))
    end

    if sub_sample_cells
        if K <= length(live_cells)
            samp_cells = sample(live_cells, K, replace = false)
            push!(samp_cell_lin_df_vec, get_counts(samp_cells, string("DT", rep, "_P", curr_P)))
        else
            push!(samp_cell_lin_df_vec, get_counts(live_cells, string("DT", rep, "_P", curr_P)))
        end
        if neutral_tracker !== nothing
            sampled = K <= length(live_cells) ? samp_cells : live_cells
            push!(neutral_tracker.sample_counts, neutral_mutation_counts(sampled, string("DT", rep, "_P", curr_P)))
        end
    end

    update_track_vec_resdmg!(kmc_out, Nvec, nS_vec, nDS_vec, nDR_vec, nR_vec, tvec, Pvec)
end

"""Phenotype-by-barcode counts for the InVivo ABM, collapsed across EG states."""
function get_pheno_counts(cells::AbstractVector{InVivoCancerCell}, rep_name::String)
    live = alive_cells(cells)
    barcodes = sort(unique(cell.barcode for cell in live))
    out = DataFrame(bc = barcodes)
    for label in ("S", "R", "E")
        out[!, Symbol(rep_name, "_", label)] = [count(cell -> cell.barcode == bc &&
            (label == "E" ? cell.E : label == "R" ? cell.R && !cell.E : !cell.R && !cell.E), live)
            for bc in barcodes]
    end
    return out
end

function _invivo_abm_sol_df(out::ResPopInVivoGrowOut; cond::String, rep::Int64)
    return DataFrame(t = out.tvec,
        nS_EG0 = out.S_EG0_vec, nS_EG1 = out.S_EG1_vec,
        nR_EG0 = out.R_EG0_vec, nR_EG1 = out.R_EG1_vec,
        nE_EG0 = out.E_EG0_vec, nE_EG1 = out.E_EG1_vec,
        nS = out.Svec, nR = out.Rvec, nE = out.Evec,
        n_EG0 = out.S_EG0_vec+out.R_EG0_vec+out.E_EG0_vec,
        n_EG1 = out.S_EG1_vec+out.R_EG1_vec+out.E_EG1_vec,
        N = out.Nvec, cond = cond, rep = rep, passage = out.Pvec)
end

function _record_invivo_abm_readout!(cells, name, lin_dfs, pheno_dfs, sub_dfs;
    full_pheno_bc::Bool = false, sub_sample_cells::Bool = false, K::Int64 = 0)
    live = alive_cells(cells)
    push!(lin_dfs, get_counts(live, name))
    full_pheno_bc && push!(pheno_dfs, get_pheno_counts(live, name))
    if sub_sample_cells
        sampled = K <= length(live) ? sample(live, K, replace = false) : live
        push!(sub_dfs, get_counts(sampled, name))
    end
    return nothing
end
