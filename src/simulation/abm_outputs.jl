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
