"""Resistance + damage ABM with ground-truth lineage tracking enabled."""
struct ResDmg_ABM_EvBC <: ABMModel
    params::ResDmgParams
    abm::ABMParams
    function ResDmg_ABM_EvBC(params::ResDmgParams, abm::ABMParams)
        validate_model_params(params)
        new(params, abm)
    end
end

"""Construct a `ResDmg_ABM_EvBC` from parameters and ABM settings."""
ResDmg_ABM_EvBC(params::ResDmgParams; abm::ABMParams = ABMParams()) = ResDmg_ABM_EvBC(params, abm)

"""Construct a `ResDmg_ABM_EvBC` from keyword arguments forwarded to `ResDmgParams`."""
ResDmg_ABM_EvBC(; abm::ABMParams = ABMParams(), kwargs...) = ResDmg_ABM_EvBC(ResDmgParams(; kwargs...), abm)

mutable struct ResDmgCellEvBC
    barcode::Float64
    DS::Bool
    DR::Bool
    R::Bool
    alive::Bool
    id::Int64
    parent_id::Int64
    birth_time::Float64
end

mutable struct ResDmgABMEvBCState
    cells::Vector{ResDmgCellEvBC}
    next_cell_id::Int64
    lineage_records::Vector{LineageRecord}
end

make_dead_resdmg_cell_evbc() = ResDmgCellEvBC(0.0, false, false, false, false, 0, 0, -1.0)

function _resdmg_pheno_label(cell::ResDmgCellEvBC)
    if cell.DS
        return "DS"
    elseif cell.DR
        return "DR"
    elseif cell.R
        return "R"
    end
    return "S"
end

function _init_lineage_state!(cells::Vector{ResDmgCellEvBC}, t0::Float64 = 0.0)
    lineage_records = LineageRecord[]
    next_cell_id = Int64(1)
    for i in eachindex(cells)
        if cells[i].alive
            cells[i].id = next_cell_id
            cells[i].parent_id = 0
            cells[i].birth_time = t0
            push!(lineage_records, LineageRecord(next_cell_id, 0, t0, "ROOT", _resdmg_pheno_label(cells[i]), cells[i].barcode))
            next_cell_id += 1
        else
            cells[i].id = 0
            cells[i].parent_id = 0
            cells[i].birth_time = -1.0
        end
    end
    return next_cell_id, lineage_records
end

function _to_evbc_cells(cells::Vector{ResDmgCell}; t0::Float64 = 0.0)
    evbc_cells = Vector{ResDmgCellEvBC}(undef, length(cells))
    for i in eachindex(cells)
        c = cells[i]
        evbc_cells[i] = ResDmgCellEvBC(c.barcode, c.DS, c.DR, c.R, c.alive, 0, 0, c.alive ? t0 : -1.0)
    end
    next_cell_id, lineage_records = _init_lineage_state!(evbc_cells, t0)
    return evbc_cells, next_cell_id, lineage_records
end

function resdmg_birth_mutate_event_evbc!(state::ResDmgABMEvBCState,
    cell_pos::Int64, birth_pos::Int64,
    mu::Float64, sig::Float64,
    phen_counts::ResDmgPhenoCounts,
    curr_t::Float64)

    cell_arr = state.cells
    if cell_pos <= 0 || cell_pos > length(cell_arr) || birth_pos <= 0 || birth_pos > length(cell_arr)
        throw(ArgumentError("Invalid cell position or birth position"))
    end

    if cell_arr[cell_pos].DS || cell_arr[cell_pos].DR
        throw(ArgumentError("Cannot perform birth event attempt on a damaged cell."))
    end

    cell_arr[birth_pos].barcode = cell_arr[cell_pos].barcode
    cell_arr[birth_pos].DS = cell_arr[cell_pos].DS
    cell_arr[birth_pos].DR = cell_arr[cell_pos].DR
    cell_arr[birth_pos].R = cell_arr[cell_pos].R
    cell_arr[birth_pos].alive = true
    child_id = state.next_cell_id
    cell_arr[birth_pos].id = child_id
    cell_arr[birth_pos].parent_id = cell_arr[cell_pos].id
    cell_arr[birth_pos].birth_time = curr_t
    parent_pheno = _resdmg_pheno_label(cell_arr[cell_pos])

    mut_p = rand()

    if cell_arr[birth_pos].R
        if mut_p < sig
            cell_arr[birth_pos].R = false
            phen_counts.Scount += 1
        else
            phen_counts.Rcount += 1
        end
    else
        if mu > mut_p
            cell_arr[birth_pos].R = true
            phen_counts.Rcount += 1
        else
            phen_counts.Scount += 1
        end
    end

    child_pheno = _resdmg_pheno_label(cell_arr[birth_pos])
    push!(state.lineage_records, LineageRecord(child_id, cell_arr[cell_pos].id, curr_t, parent_pheno, child_pheno, cell_arr[birth_pos].barcode))
    state.next_cell_id += 1
end

function _core_grow_kill_abm!(
    state::ResDmgABMEvBCState,
    params::ResDmgParams,
    sim::ABMSimParams;
    treat::Bool = false
)
    @assert 0.0 <= sim.t_frac <= 1.0 "t_frac must be between 0 and 1."
    @assert params.drug_effect in RESDMG_DRUG_EFFECTS "drug effect can only take :b, :d, or :c."
    @assert sim.R_real in ["b", "d", "l"] "R_real can only take 'b', 'd' or 'l' as a value."

    cells = state.cells
    de = sim.drug_effect
    if de == :b
        @assert params.Dc <= params.b "When drug_effect = :b, Dc must be <= b."
    end

    bmax = params.b
    omemax_S = params.ome
    omemax_R = params.ome * (1 - params.psi)
    zetmax_S = params.zet_S
    zetmax_R = params.zet_R
    zetmax = max(zetmax_S, zetmax_R)
    omemax = max(omemax_S, omemax_R)
    lam = bmax - params.d

    dmax_no_treat = if params.del > 0.0 && sim.R_real == "d"
        params.d + (lam * params.del)
    else
        params.d
    end

    dmax_treat = dmax_no_treat + params.Dc

    drug_concs = nothing
    if treat
        drug_concs = drug_treat_concs(sim.tmax, params.k, params.Dc,
                                      sim.treat_ons, sim.treat_offs, sim.dt_save_at)
        bdmax = bmax + dmax_treat
    else
        bdmax = bmax + dmax_no_treat
    end
    bdzetmax = bdmax + zetmax
    bdzetomemax = bdmax + omemax + zetmax

    t = sim.t0
    t_rec_change = sim.tmax * sim.t_frac
    t_rec = t_rec_change
    tvec = Float64[t]

    order_cells!(cells, make_dead_resdmg_cell_evbc)
    live_vec = live_positions(cells)
    dead_vec = dead_positions(cells)

    Nt = n_alive(cells)
    DScount = Int64(sum(cell -> cell.DS && cell.alive, cells))
    DRcount = Int64(sum(cell -> cell.DR && cell.alive, cells))
    Rcount = Int64(sum(cell -> cell.R && cell.alive, cells))
    Scount = Int64(sum(cell -> !cell.DS && !cell.DR && !cell.R && cell.alive, cells))

    phen_counts = ResDmgPhenoCounts(Scount, DScount, DRcount, Rcount)

    Nvec = Int64[Nt]
    Svec = Int64[phen_counts.Scount]
    DSvec = Int64[phen_counts.DScount]
    DRvec = Int64[phen_counts.DRcount]
    Rvec = Int64[phen_counts.Rcount]
    Pvec = Int64[sim.Passage]

    while t <= sim.tmax
        Nt > 0 || error("ABM core loop encountered non-positive live population (Nt=$(Nt)).")
        any(!c.alive for c in cells) || error("ABM core loop has no available dead slots for births; increase Nbuff.")

        live_pos = rand(1:length(cells))
        dead_pos = rand(1:length(cells))

        while live_vec[live_pos] == 0
            live_pos = rand(1:length(cells))
        end

        while dead_vec[dead_pos] == 0
            dead_pos = rand(1:length(cells))
        end

        curr_rconc = 0.0
        if treat
            curr_rconc = curr_rc(t, drug_concs)
            dmax_window = curr_rconc > 0.0 ? dmax_treat : dmax_no_treat
            ran_max = curr_rconc > 0.0 ? bdzetomemax : bdzetmax
            ran = rand(Uniform(0, ran_max))
            denom = ran_max * Nt
            denom > 0.0 || error("ABM core loop encountered non-positive time-step denominator (rate*Nt=$(denom)).")
            dt = -log(rand()) / denom
            t += dt
        else
            dmax_window = dmax_no_treat
            ran = rand(Uniform(0, bdzetmax))
            denom = bdzetmax * Nt
            denom > 0.0 || error("ABM core loop encountered non-positive time-step denominator (rate*Nt=$(denom)).")
            dt = -log(rand()) / denom
            t += dt
        end

        if t >= t_rec
            push!(Nvec, Nt)
            push!(tvec, t)
            push!(Svec, phen_counts.Scount)
            push!(DSvec, phen_counts.DScount)
            push!(DRvec, phen_counts.DRcount)
            push!(Rvec, phen_counts.Rcount)
            push!(Pvec, sim.Passage)
            t_rec += t_rec_change
        end

        if t >= sim.tmax
            t = sim.tmax
            break
        end

        if Nt >= sim.Nmax
            push!(Nvec, Nt)
            push!(tvec, t)
            push!(Svec, phen_counts.Scount)
            push!(DSvec, phen_counts.DScount)
            push!(DRvec, phen_counts.DRcount)
            push!(Rvec, phen_counts.Rcount)
            push!(Pvec, sim.Passage)
            break
        end

        if params.del > 0.0
            if cells[live_pos].R
                cell_b, cell_d = if sim.R_real == "b"
                    (params.b - (lam * params.del), params.d)
                elseif sim.R_real == "d"
                    (params.b, params.d + (lam * params.del))
                else
                    (params.b * (1 - params.del), params.d * (1 - params.del))
                end
            else
                cell_b = params.b
                cell_d = params.d
            end
            cell_b = max(cell_b, 0.0)
            cell_d = max(cell_d, 0.0)
        else
            cell_b = params.b
            cell_d = params.d
        end

        if cells[live_pos].DS || cells[live_pos].DR
            cell_b = 0.0
            cell_d = params.d + params.Dc
        elseif treat
            curr_dconc = curr_dc(t, drug_concs)
            if de == :d
                cell_d += if cells[live_pos].R
                    curr_dconc * (1 - params.psi)
                else
                    curr_dconc
                end
            elseif de == :b
                cell_b -= if cells[live_pos].R
                    curr_dconc * (1 - params.psi)
                else
                    curr_dconc
                end
            elseif de == :c
                cell_b -= if cells[live_pos].R
                    curr_dconc * (1 - params.psi)
                else
                    curr_dconc
                end
                if cell_b < 0.0
                    negative_surplus = abs(cell_b)
                    cell_b = 0.0
                    cell_d += negative_surplus
                end
            end
        end

        cell_b *= (1 - (Nt / sim.Cc))
        cell_d *= (1 - (Nt / sim.Cc))

        if treat
            if cells[live_pos].R
                cell_ome = params.ome * (1 - params.psi) * curr_rconc
            elseif cells[live_pos].DS || cells[live_pos].DR
                cell_ome = 0.0
            else
                cell_ome = params.ome * curr_rconc
            end
        else
            cell_ome = 0.0
        end

        cell_ome <= omemax || error("Calculated cell_ome value $(cell_ome) exceeds maximum possible value $(omemax). Check parameter values and drug effect settings.")

        if cells[live_pos].DS
            cell_zet = params.zet_S
        elseif cells[live_pos].DR
            cell_zet = params.zet_R
        else
            cell_zet = 0.0
        end

        cell_zet <= zetmax || error("Calculated cell_zet value $(cell_zet) exceeds maximum possible value $(zetmax). Check parameter values.")

        if ran < cell_b
            resdmg_birth_mutate_event_evbc!(state, live_pos, dead_pos,
                                            params.mu, params.sig,
                                            phen_counts, t)

            live_vec[dead_pos] = 1
            dead_vec[dead_pos] = 0
            Nt += 1
        end

        if bmax <= ran < bmax + cell_d
            resdmg_death_event!(cells, live_pos, phen_counts)
            live_vec[live_pos] = 0
            dead_vec[live_pos] = 1
            Nt -= 1
        end

        if (bmax + dmax_window) <= ran < (bmax + dmax_window + cell_zet)
            resdmg_repair_event!(cells, live_pos, phen_counts)
        end

        if treat && ((bmax + dmax_window + zetmax) <= ran < (bmax + dmax_window + zetmax + cell_ome))
            resdmg_damage_event!(cells, live_pos, phen_counts)
        end

        if Nt <= 0
            push!(Nvec, Nt)
            push!(tvec, t)
            push!(Svec, phen_counts.Scount)
            push!(DSvec, phen_counts.DScount)
            push!(DRvec, phen_counts.DRcount)
            push!(Rvec, phen_counts.Rcount)
            push!(Pvec, sim.Passage)
            break
        end
    end

    fin_t = round(t, digits = 4)
    return ResDmgGrowOut(Nvec, tvec, Svec, DSvec, DRvec, Rvec, Pvec, fin_t)
end

function run_model_core_abm(model::ResDmg_ABM_EvBC, state::ResDmgABMEvBCState, sim::ABMSimParams; treat::Bool = false)
    return _core_grow_kill_abm!(state, model.params, sim; treat = treat)
end
