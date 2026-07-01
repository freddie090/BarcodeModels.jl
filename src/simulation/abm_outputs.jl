"""Append ResPop ABM trajectory and barcode outputs for one passage."""
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

"""Append in vivo EG-stratified ABM counts to output vectors."""
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

"""Append an in vivo ABM segment to accumulated solution vectors."""
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

"""Build the in vivo ABM population trajectory DataFrame."""
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

"""Append ResDmg ABM trajectory and barcode outputs for one passage."""
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
