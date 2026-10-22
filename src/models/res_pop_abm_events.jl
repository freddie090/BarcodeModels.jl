# Standard ABM birth and shared events

"""Apply a ResPop ABM birth with a possible daughter phenotype change."""
function birth_mutate_event!(cell_arr,
    cell_pos::Int64, birth_pos::Int64,
    mu::Float64, sig::Float64, al::Float64,
    phen_counts)

    if cell_pos <= 0 || cell_pos > length(cell_arr) || birth_pos <= 0 || birth_pos > length(cell_arr)
        throw(ArgumentError("Invalid cell position or birth position"))
    end

    cell_arr[birth_pos].barcode = cell_arr[cell_pos].barcode
    cell_arr[birth_pos].R = cell_arr[cell_pos].R
    cell_arr[birth_pos].E = cell_arr[cell_pos].E
    cell_arr[birth_pos].alive = true

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

# Independent conversion changes phenotype only, retaining barcode and cell identity.
# AbstractVector also serves the existing EvBC cell workflow.
function respop_S_to_R_transition_event!(cells::AbstractVector, cell_pos::Int64,
    phen_counts::PhenoCounts)
    cell = cells[cell_pos]
    cell.alive && !cell.R && !cell.E || return
    cell.R = true
    phen_counts.Scount -= 1
    phen_counts.Rcount += 1
end

function respop_R_to_S_transition_event!(cells::AbstractVector, cell_pos::Int64,
    phen_counts::PhenoCounts)
    cell = cells[cell_pos]
    cell.alive && cell.R && !cell.E || return
    cell.R = false
    phen_counts.Rcount -= 1
    phen_counts.Scount += 1
end

function respop_R_to_E_transition_event!(cells::AbstractVector, cell_pos::Int64,
    phen_counts::PhenoCounts)
    cell = cells[cell_pos]
    cell.alive && cell.R && !cell.E || return
    cell.R = false
    cell.E = true
    phen_counts.Rcount -= 1
    phen_counts.Ecount += 1
end

"""Mark a ResPop-family ABM cell dead and decrement its phenotype count."""
function death_event!(cell_arr, cell_pos::Int64, phen_counts)
    if cell_pos <= 0 || cell_pos > length(cell_arr)
        throw(ArgumentError("Invalid cell position"))
    end

    cell_arr[cell_pos].alive = false

    if cell_arr[cell_pos].E
        phen_counts.Ecount -= 1
    elseif cell_arr[cell_pos].R
        phen_counts.Rcount -= 1
    else
        phen_counts.Scount -= 1
    end
end

# Instantaneous treatment events (standard ABM and EvBC)

"""Apply one instantaneous kill, preserving survivor metadata and lineage records."""
function _instant_kill_event!(cells, phen_counts::PhenoCounts, live_vec, dead_vec,
    event::InstantKillState)
    event.applied && error("Instantaneous killing event has already been applied.")
    abm = event.settings
    if abm.inst_kill_mode == "neutral_mutation"
        event.mutation_id, event.carrier_fraction = _instant_resistance_mutation(cells, abm)
    end
    event.N_before = n_alive(cells)
    killed = 0
    for i in eachindex(cells)
        cell = cells[i]
        cell.alive || continue
        resistant = abm.inst_kill_mode == "phenotype" ? cell.R || cell.E :
            abm.inst_kill_mode == "neutral_mutation" ? event.mutation_id in cell.muts : false
        p = abm.inst_drug_kill_prob * (resistant ? abm.inst_R_strength : 1.0)
        if p == 1.0 || (p > 0.0 && rand() < p)
            death_event!(cells, i, phen_counts)
            live_vec[i] = 0
            dead_vec[i] = 1
            killed += 1
        end
    end
    event.N_killed = killed
    event.N_after = event.N_before - killed
    event.applied = true
    return event.N_after
end

"""Record both sides of the instantaneous discontinuity at its exact time."""
function _record_instant_kill!(cells, phen_counts, live_vec, dead_vec, event,
    sim, Nvec, tvec, Svec, Rvec, Evec, Pvec)
    time = event.settings.inst_kill_time
    push!(Nvec, n_alive(cells))
    push!(tvec, time)
    push!(Svec, phen_counts.Scount)
    push!(Rvec, phen_counts.Rcount)
    push!(Evec, phen_counts.Ecount)
    push!(Pvec, sim.Passage)
    Nt = _instant_kill_event!(cells, phen_counts, live_vec, dead_vec, event)
    push!(Nvec, Nt)
    push!(tvec, time)
    push!(Svec, phen_counts.Scount)
    push!(Rvec, phen_counts.Rcount)
    push!(Evec, phen_counts.Ecount)
    push!(Pvec, sim.Passage)
    return Nt
end

# EvBC birth event

function _birth_mutate_event_evbc!(state::ResPopABMEvBCState,
    cell_pos::Int64, birth_pos::Int64,
    mu::Float64, sig::Float64, al::Float64,
    phen_counts::PhenoCounts,
    curr_t::Float64)

    cell_arr = state.cells
    if cell_pos <= 0 || cell_pos > length(cell_arr) || birth_pos <= 0 || birth_pos > length(cell_arr)
        throw(ArgumentError("Invalid cell position or birth position"))
    end

    mother_id = cell_arr[cell_pos].id
    mother_barcode = cell_arr[cell_pos].barcode
    mother_R = cell_arr[cell_pos].R
    mother_E = cell_arr[cell_pos].E
    mother_pheno = _respop_pheno_label(cell_arr[cell_pos])

    daughter_a_id = state.next_cell_id
    cell_arr[cell_pos].id = daughter_a_id
    cell_arr[cell_pos].parent_id = mother_id
    cell_arr[cell_pos].birth_time = curr_t
    push!(state.lineage_records, LineageRecord(daughter_a_id, mother_id, curr_t, mother_pheno, mother_pheno, mother_barcode))
    state.next_cell_id += 1

    cell_arr[birth_pos].barcode = mother_barcode
    cell_arr[birth_pos].R = mother_R
    cell_arr[birth_pos].E = mother_E
    cell_arr[birth_pos].alive = true
    daughter_b_id = state.next_cell_id
    cell_arr[birth_pos].id = daughter_b_id
    cell_arr[birth_pos].parent_id = mother_id
    cell_arr[birth_pos].birth_time = curr_t

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

    child_pheno = _respop_pheno_label(cell_arr[birth_pos])
    push!(state.lineage_records, LineageRecord(daughter_b_id, mother_id, curr_t, mother_pheno, child_pheno, cell_arr[birth_pos].barcode))
    state.next_cell_id += 1
end

