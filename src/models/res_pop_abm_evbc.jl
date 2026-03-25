"""Lineage record for evolving-barcode simulations."""
struct LineageRecord
    id::Int64
    parent_id::Int64
    birth_time::Float64
end

"""Resistance population ABM with ground-truth lineage tracking enabled."""
struct ResPop_ABM_EvBC <: ABMModel
    params::ResPopParams
    abm::ABMParams
    function ResPop_ABM_EvBC(params::ResPopParams, abm::ABMParams)
        validate_model_params(params)
        new(params, abm)
    end
end

"""Construct a `ResPop_ABM_EvBC` from parameters and ABM settings."""
ResPop_ABM_EvBC(params::ResPopParams; abm::ABMParams = ABMParams()) = ResPop_ABM_EvBC(params, abm)

"""Construct a `ResPop_ABM_EvBC` from keyword arguments forwarded to `ResPopParams`."""
ResPop_ABM_EvBC(; abm::ABMParams = ABMParams(), kwargs...) = ResPop_ABM_EvBC(ResPopParams(; kwargs...), abm)

mutable struct CancerCellEvBC
    barcode::Float64
    R::Bool
    E::Bool
    alive::Bool
    id::Int64
    parent_id::Int64
    birth_time::Float64
end

mutable struct ResPopABMEvBCState
    cells::Vector{CancerCellEvBC}
    next_cell_id::Int64
    lineage_records::Vector{LineageRecord}
end

make_dead_cell_evbc() = CancerCellEvBC(0.0, false, false, false, 0, 0, -1.0)

function _init_lineage_state!(cells::Vector{CancerCellEvBC}, t0::Float64 = 0.0)
    lineage_records = LineageRecord[]
    next_cell_id = Int64(1)
    for i in eachindex(cells)
        if cells[i].alive
            cells[i].id = next_cell_id
            cells[i].parent_id = 0
            cells[i].birth_time = t0
            push!(lineage_records, LineageRecord(next_cell_id, 0, t0))
            next_cell_id += 1
        else
            cells[i].id = 0
            cells[i].parent_id = 0
            cells[i].birth_time = -1.0
        end
    end
    return next_cell_id, lineage_records
end

function _to_evbc_cells(cells::Vector{CancerCell}; t0::Float64 = 0.0)
    evbc_cells = Vector{CancerCellEvBC}(undef, length(cells))
    for i in eachindex(cells)
        c = cells[i]
        evbc_cells[i] = CancerCellEvBC(c.barcode, c.R, c.E, c.alive, 0, 0, c.alive ? t0 : -1.0)
    end
    next_cell_id, lineage_records = _init_lineage_state!(evbc_cells, t0)
    return evbc_cells, next_cell_id, lineage_records
end

function _birth_mutate_event_evbc!(state::ResPopABMEvBCState,
    cell_pos::Int64, birth_pos::Int64,
    mu::Float64, sig::Float64, al::Float64,
    phen_counts::PhenoCounts,
    curr_t::Float64)

    cell_arr = state.cells
    if cell_pos <= 0 || cell_pos > length(cell_arr) || birth_pos <= 0 || birth_pos > length(cell_arr)
        throw(ArgumentError("Invalid cell position or birth position"))
    end

    cell_arr[birth_pos].barcode = cell_arr[cell_pos].barcode
    cell_arr[birth_pos].R = cell_arr[cell_pos].R
    cell_arr[birth_pos].E = cell_arr[cell_pos].E
    cell_arr[birth_pos].alive = true
    cell_arr[birth_pos].id = state.next_cell_id
    cell_arr[birth_pos].parent_id = cell_arr[cell_pos].id
    cell_arr[birth_pos].birth_time = curr_t
    push!(state.lineage_records, LineageRecord(state.next_cell_id, cell_arr[cell_pos].id, curr_t))
    state.next_cell_id += 1

    mut_p = rand()

    if cell_arr[birth_pos].E
        phen_counts.Ecount += 1
    elseif cell_arr[birth_pos].R
        if mut_p < sig
            cell_arr[birth_pos].R = false
            phen_counts.Scount += 1
        elseif sig <= mut_p < (sig + al)
            cell_arr[birth_pos].R = false
            cell_arr[birth_pos].E = true
            phen_counts.Ecount += 1
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
end

function _core_grow_kill_abm!(
    state::ResPopABMEvBCState,
    params::ResPopParams,
    sim::ABMSimParams;
    treat::Bool = false
)
    @assert 0.0 <= sim.t_frac <= 1.0 "t_frac must be between 0 and 1."
    @assert params.drug_effect in RESPOP_DRUG_EFFECTS "drug effect can only take :b, :d, or :c."
    @assert sim.R_real in ["b", "d", "l"] "R_real can only take 'b', 'd' or 'l' as a value."

    cells = state.cells
    de = sim.drug_effect
    if de == :b
        @assert params.Dc <= params.b "When drug_effect = :b, Dc must be <= b."
    end

    bmax = params.b
    dmax = params.d
    lam = bmax - dmax

    if params.del > 0.0
        if sim.R_real == "d"
            dmax = dmax + (lam * params.del)
            c_bdmax = bmax + dmax
        else
            c_bdmax = bmax + dmax
        end
    else
        c_bdmax = bmax + dmax
    end

    drug_concs = nothing
    if treat
        drug_concs = drug_treat_concs(sim.tmax, params.k, params.Dc,
                                      sim.treat_ons, sim.treat_offs, sim.dt_save_at)
        if de == :d || de == :c
            dmax = params.d + maximum(drug_concs["dconc"])
            t_bdmax = bmax + dmax
        else
            t_bdmax = bmax + dmax
        end
    else
        t_bdmax = bmax + dmax
    end

    bdmax = maximum([c_bdmax, t_bdmax])

    t = sim.t0
    t_rec_change = sim.tmax * sim.t_frac
    t_rec = t_rec_change
    tvec = Float64[t]

    order_cells!(cells, make_dead_cell_evbc)
    live_vec = live_positions(cells)
    dead_vec = dead_positions(cells)

    Nt = n_alive(cells)
    Rcount = Int64(sum(cell -> cell.R && cell.alive, cells))
    Ecount = Int64(sum(cell -> cell.E && cell.alive, cells))
    Scount = Int64(sum(cell -> !cell.R && !cell.E && cell.alive, cells))
    phen_counts = PhenoCounts(Scount, Rcount, Ecount)

    Nvec = Int64[Nt]
    Svec = Int64[phen_counts.Scount]
    Rvec = Int64[phen_counts.Rcount]
    Evec = Int64[phen_counts.Ecount]
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

        ran = rand(Uniform(0, bdmax))
        denom = bdmax * Nt
        denom > 0.0 || error("ABM core loop encountered non-positive time-step denominator (bdmax*Nt=$(denom)).")
        dt = -log(rand()) / denom
        t += dt

        if t >= t_rec
            push!(Nvec, Nt)
            push!(tvec, t)
            push!(Svec, phen_counts.Scount)
            push!(Rvec, phen_counts.Rcount)
            push!(Evec, phen_counts.Ecount)
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
            push!(Rvec, phen_counts.Rcount)
            push!(Evec, phen_counts.Ecount)
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
            cell_b = max(cell_b, 0)
            cell_d = max(cell_d, 0)
        else
            cell_b = params.b
            cell_d = params.d
        end

        if treat
            curr_dconc = curr_dc(t, drug_concs)
            if de == :d
                cell_d += if cells[live_pos].R || cells[live_pos].E
                    curr_dconc * (1 - params.psi)
                else
                    curr_dconc
                end
            elseif de == :b
                cell_b -= if cells[live_pos].R || cells[live_pos].E
                    curr_dconc * (1 - params.psi)
                else
                    curr_dconc
                end
            elseif de == :c
                cell_b -= if cells[live_pos].R || cells[live_pos].E
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

        if ran < cell_b
            if treat
                curr_rconc = curr_rc(t, drug_concs)
                al_scal = params.al * curr_rconc
                _birth_mutate_event_evbc!(state, live_pos, dead_pos,
                                          params.mu, params.sig, al_scal,
                                          phen_counts, t)
            else
                _birth_mutate_event_evbc!(state, live_pos, dead_pos,
                                          params.mu, params.sig, 0.0,
                                          phen_counts, t)
            end

            live_vec[dead_pos] = 1
            dead_vec[dead_pos] = 0
            Nt += 1
        end

        if bmax <= ran < bmax + cell_d
            death_event!(cells, live_pos, phen_counts)
            live_vec[live_pos] = 0
            dead_vec[live_pos] = 1
            Nt -= 1
        end

        if Nt <= 0
            push!(Nvec, Nt)
            push!(tvec, t)
            push!(Svec, phen_counts.Scount)
            push!(Rvec, phen_counts.Rcount)
            push!(Evec, phen_counts.Ecount)
            push!(Pvec, sim.Passage)
            break
        end
    end

    fin_t = round(t, digits = 4)
    return GrowOut(Nvec, tvec, Svec, Rvec, Evec, Pvec, fin_t)
end

function run_model_core_abm(model::ResPop_ABM_EvBC, state::ResPopABMEvBCState, sim::ABMSimParams; treat::Bool = false)
    return _core_grow_kill_abm!(state, model.params, sim; treat = treat)
end
