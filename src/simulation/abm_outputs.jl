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

"""Record ResDmg ABM trajectory, barcode and optional subsample outputs."""
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

